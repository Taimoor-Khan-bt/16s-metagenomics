# Config Inheritance and Merging Utility
# Allows project configs to inherit from base_config.yaml

merge_configs <- function(base_path, override_path) {
  # Load base config
  if (!file.exists(base_path)) {
    stop("Base config not found: ", base_path)
  }
  base <- yaml::read_yaml(base_path)
  
  # Load override config
  if (!file.exists(override_path)) {
    stop("Override config not found: ", override_path)
  }
  override <- yaml::read_yaml(override_path)
  
  # Deep merge function
  deep_merge <- function(base_list, override_list) {
    if (!is.list(base_list) || !is.list(override_list)) {
      # If not both lists, override takes precedence
      return(override_list)
    }
    
    # Start with base
    result <- base_list
    
    # Merge each element from override
    for (name in names(override_list)) {
      if (name %in% names(result) && 
          is.list(result[[name]]) && 
          is.list(override_list[[name]])) {
        # Recursively merge nested lists
        result[[name]] <- deep_merge(result[[name]], override_list[[name]])
      } else {
        # Direct override
        result[[name]] <- override_list[[name]]
      }
    }
    
    return(result)
  }
  
  # Perform merge
  merged <- deep_merge(base, override)
  
  return(merged)
}

# Check if config has inheritance directive
has_inheritance <- function(config_path) {
  if (!file.exists(config_path)) {
    return(FALSE)
  }
  
  config <- yaml::read_yaml(config_path)
  return(!is.null(config$inherit_from))
}

# Load config with optional inheritance
load_config_with_inheritance <- function(config_path, base_config_dir = "config") {
  config <- yaml::read_yaml(config_path)
  
  # Check for inheritance directive
  if (!is.null(config$inherit_from)) {
    inherit_from <- config$inherit_from
    
    # Resolve base config path
    if (!grepl("^/", inherit_from)) {
      # Relative path - resolve from config directory
      base_path <- file.path(base_config_dir, inherit_from)
    } else {
      base_path <- inherit_from
    }
    
    message("[config] Inheriting from: ", base_path)
    
    # Merge configs
    merged <- merge_configs(base_path, config_path)
    
    # Remove inheritance directive from final config
    merged$inherit_from <- NULL
    
    return(merged)
  }
  
  # No inheritance - return as-is
  return(config)
}

# Validate that required fields are present after inheritance
validate_inherited_config <- function(config, required_sections) {
  missing <- character()
  
  for (section in required_sections) {
    if (is.null(config[[section]])) {
      missing <- c(missing, section)
    }
  }
  
  if (length(missing) > 0) {
    stop("Config missing required sections after inheritance: ", 
         paste(missing, collapse = ", "))
  }
  
  return(invisible(TRUE))
}

# Example usage:
# cfg <- load_config_with_inheritance("config/my_project.yaml")
# validate_inherited_config(cfg, c("project", "io", "amplicon"))
